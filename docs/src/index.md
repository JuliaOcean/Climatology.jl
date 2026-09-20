
# Climatology.jl

Climatology.jl provides analysis tools and workflows for accessing and deriving 
climatologies from impactful ocean and climate data sets. It serves as a shared 
data-access and diagnostics layer: gridded fields are downloaded lazily, read into 
standard Julia array/table structures, and processed into derived quantities, maps, 
and time series — reproducibly, from interactive Pluto notebooks or scripted workflows.

_It is in early development stage; breaking changes remain likely._

## Supported Data Sets

| Data set | Description | Status |
|---|---|---|
| [OISST](@ref Physical-Oceanography) | NOAA sea surface temperature, incl. anomalies & marine heat wave detection | most developed |
| [ECCO](@ref Physical-Oceanography) | ocean state estimate: transports, climate indices, sub-surface temperature, etc. | most developed |
| OCCA (1 & 2) | ocean state estimates or climatologies | supported |
| [CBIOMES](@ref Marine-Ecosystems) | ocean color & biogeochemistry from Darwin3 | supported |
| MITprof | in-situ temperature/salinity profile climatology | supported |
| [HadIOD](@ref Other-Notebooks) | in-situ T/S observational database | supported |
| Satellite altimetry / SSH | sea level anomaly maps (NASA/PODAAC, CMEMS) | supported |
| NSLCT | NASA sea level time series & maps | supported |

## Design

- **Downloads**: data sets are fetched lazily, on demand, via 
  [Scratch.jl](https://github.com/JuliaPackaging/Scratch.jl), 
  [DataDeps.jl](https://github.com/oxinabox/DataDeps.jl), and 
  [Dataverse.jl](https://gdcc.github.io/Dataverse.jl/stable/). 
  See [API](@ref) for artifact paths and download functions.
- **Plotting**: plot recipes are provided via a `Makie.jl` package extension — 
  Makie is an optional, not a hard, dependency.
- **Diagnostics & options**: each plot/analysis type is driven by a typed, 
  extensible set of options that parametrizes both the underlying computation 
  and the resulting plot.
- **Time series analysis**: trends and related statistics are computed via 
  [GLM.jl](https://github.com/JuliaStats/GLM.jl).
- **Grids**: non-regular model/observational grids (e.g. LLC90) are handled via 
  [MeshArrays.jl](https://github.com/JuliaClimate/MeshArrays.jl).
- **Derived products**: workflows write intermediate/derived data products back 
  to disk, including NetCDF, for reuse without recomputing from raw sources.
- **Notebooks**: most examples are distributed as interactive 
  [Pluto.jl](https://github.com/fonsp/Pluto.jl) notebooks (see [Examples](@ref)).

## Related Packages

- [MeshArrays.jl](https://github.com/JuliaClimate/MeshArrays.jl) and [Drifters.jl](https://github.com/JuliaClimate/Drifters.jl) depend on 
  Climatology.jl for its example suite and test cases.
- [ArgoData.jl](https://github.com/JuliaOcean/ArgoData.jl) consumes Climatology.jl, 
  via a package extension, to download gridded climatologies used for sampling 
  MITprof profiles. ArgoData.jl also independently derives gridded climatologies 
  from irregularly-distributed Argo profile data via geospatial statistics; some 
  of that workflow may be consolidated into Climatology.jl in a future refactor.

## Scope Note

Climatology.jl serves and analyzes *already-computed* reanalysis/observational 
output. It does not run ocean/climate models or perform adjoint/inverse modeling — 
for that, see [MITgcm.jl](https://gaelforget.github.io/MITgcm.jl/dev/) and 
[ECCO.jl](https://github.com/gaelforget/ECCO.jl).

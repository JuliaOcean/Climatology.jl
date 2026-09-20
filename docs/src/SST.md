# SST

The workflow presented here parallels the [ECCO](@ref) workflow, applied to 
NOAA's OISST sea surface temperature product.

- set up for running analyses of `SST` (OISST) data.
- run a diagnostic (e.g. `by_time`, `MHW`, `local_and_global`, `TimeLat`) over the 
  configured time period.

## Supported `plot_type`s

<!-- TODO: confirm this list against src/types.jl's SST section -->
| `plot_type` | Description |
|---|---|
| `:by_time` | time series over a chosen `period` |
| `:MHW` | marine heat wave detection |
| `:local_and_global` | local vs. global mean comparison |
| `:TimeLat` | time-latitude Hovmöller-style section |
| `:by_year` | year-by-year comparison |
| `:map_base` / `:map` | map view |

## I/O and Diagnostics Reference

```@autodocs
Modules = [Climatology.SST_FILES,SST_processing,SST_coarse_grain,SST_timeseries,SST_scenarios]
```

## Plotting (Makie extension)

```@autodocs
Modules = [ClimatologyMakieExt.SST_plots,ClimatologyMakieExt.ERA5_plot]
```

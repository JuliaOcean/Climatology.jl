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
| `:ersst_map` | ERSST monthly anomaly map (data from `SST_FILES.ersst_map_data`) |
| `:map_point` | daily map with a 10° box marker, plus the box series beside it (`map_data`, `point`, `timeseries`, `lon`, `lat`, `show_anom`) |

## I/O and Diagnostics Reference

```@autodocs
Modules = [Climatology.SST_FILES,SST_processing,SST_coarse_grain,SST_timeseries,SST_scenarios]
```

## Plotting (Makie extension)

```@autodocs
Modules = [ClimatologyMakieExt.SST_plots,ClimatologyMakieExt.ERA5_plot]
```

## Surface Fluxes (ERA5)

```@autodocs
Modules = [Climatology.ERA5]
```


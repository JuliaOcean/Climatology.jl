# Guide

This page describes the mental model behind Climatology.jl's design — how data 
gets in, how it's organized once in memory, and how computation connects to 
plotting. For function-level reference, see [API](@ref). For the ECCO-specific 
workflow, see [ECCO](@ref).

## Data flow

1. **Download** — a data set is fetched lazily, on first use, into a package-managed 
   scratch space or DataDep folder (see [API](@ref) for the mechanics).
2. **Read** — files (mostly NetCDF) are read into standard Julia array/table 
   structures, or into `MeshArrays.jl` structures when the underlying grid is 
   non-regular (e.g. the LLC90 grid used by ECCO).
3. **Diagnose** — a diagnostic/plot type is selected and configured via an 
   *options* object (see below), which drives both a computation and, optionally, 
   a plot.
4. **Derive & persist** — results can be written back to disk (often NetCDF) as 
   intermediate/derived products, so downstream steps don't need to recompute 
   from raw sources.

## The options pattern

Each supported diagnostic or plot type (e.g. an ECCO overturning circulation 
plot, an SST marine-heat-wave map) is identified by a `plot_type` symbol and 
configured through a `NamedTuple` of options. Three small functions manage this:

- `default_options(::Val{plot_type})` returns the baseline options for a given 
  `plot_type`, dispatched on `Val` so each plot type can declare its own defaults 
  (e.g. a default `period`, or no default at all if a field is required at 
  construction time).
- `getopt`/`setopt` read and update fields on an options `NamedTuple` defensively — 
  falling back gracefully instead of raising a `KeyError` when an optional field 
  is absent (e.g. `year_range`, used to derive an averaging window or axis range 
  from `period` when `years_to_display` hasn't been set).

The same options object flows through both the computation function (e.g. 
averaging a field over a chosen time window) and, if `Makie.jl` is loaded, the 
corresponding plot recipe — so a single configuration governs what gets computed 
*and* how it's displayed.

This pattern is shared across `ECCOdiag` and `SSTdiag` plot types; see 
[ECCO](@ref) and [SST](@ref) for the concrete list of supported `plot_type`s and 
their options.

## Plotting via package extension

Plotting recipes live in a `Makie.jl` package extension, not in Climatology.jl's 
core code. This keeps `Makie.jl` — a fairly heavy dependency — optional: 
Climatology.jl can be used purely for data access and computation without ever 
loading a plotting backend. Recipes become available automatically once the user 
does `using Makie` (or a Makie backend such as `CairoMakie`) alongside 
Climatology.jl.

## Grids

Regular lon-lat grids are handled with plain Julia arrays. Non-regular 
model/observational grids — such as the LLC90 grid used by ECCO — are handled via 
[MeshArrays.jl](https://github.com/JuliaClimate/MeshArrays.jl), which provides 
domain decomposition and the array-of-arrays abstraction needed to treat a 
multi-facet grid as a single object.

## Time series analysis

Trends and related statistics (e.g. linear trend, confidence intervals) are 
computed via [GLM.jl](https://github.com/JuliaStats/GLM.jl), applied consistently 
across the time-series-producing diagnostics (e.g. global means, regional means).

## Where this connects to other packages

Climatology.jl is a shared data-access and diagnostics layer. `ArgoData.jl` 
builds on it (via a package extension) to obtain gridded climatologies for 
profile sampling; `MeshArrays.jl` and `Drifters.jl` depend on it for their 
example suites and test cases. See [index](@ref) for the full picture.


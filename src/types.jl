
##

abstract type AbstractClimateDiagnostic <: Any end

##

Base.@kwdef struct ECCOdiag <: AbstractClimateDiagnostic
    path :: String = tempdir()
    name :: String = "unknown"
    options :: NamedTuple = NamedTuple()
    data :: AbstractArray = []
end

###

default_options(plot_type::Symbol) = default_options(Val(plot_type))
default_options(::Val{T}) where T = (plot_type=T,)

## SSTdiag
default_options(::Val{:by_year})          = (plot_type=:by_year,)
default_options(::Val{:by_time})          = (plot_type=:by_time, period=(1982,2024), show_anom=true, show_clim=true)
default_options(::Val{:MHW})              = (plot_type=:MHW, period=(1982,2024))
default_options(::Val{:local_and_global}) = (plot_type=:local_and_global, period=(1982,2024), ylims=(-2.5,2.5))
default_options(::Val{:map_base})         = (plot_type=:map_base,)
default_options(::Val{:map})              = (plot_type=:map,)
default_options(::Val{:TimeLat})          = (plot_type=:TimeLat, period=(1982,2024), ylims=(-90,90), clip_to_range=true)

"""
    default_options(::Val{T}) where T

Return the baseline options `NamedTuple` for a given `plot_type` symbol `T`,
dispatched via `Val` so that each plot type can declare its own defaults.

The generic fallback returns just `(plot_type=T,)`. Specific plot types
override this method (dispatching on their own `Val{:plot_type_name}`) to
add fields such as `period` or `years_to_display` — or to omit them
entirely when no sensible default exists (e.g. required opaque payloads
like a grid `Γ`, supplied only at construction time).

# Examples
```julia
default_options(Val(:ECCO_OHT1))
# (plot_type = :ECCO_OHT1, period = (1992, 2011), years_to_display = nothing)
```

See also [`getopt`](@ref), [`setopt`](@ref), [`year_range`](@ref).
"""
function default_options end

default_options(::Val{:ECCO_TimeLat}) = (
    plot_type=:ECCO_TimeLat, select_method=0, period=(1992,2011),
    level=1, ylims=(-90,90), colormap_factor=1, years_to_display=nothing,
)
default_options(::Val{:ECCO_TimeLatAnom}) =
    merge(default_options(Val(:ECCO_TimeLat)), (plot_type=:ECCO_TimeLatAnom, select_method=1))
default_options(::Val{:ECCO_GlobalMean}) = (plot_type=:ECCO_GlobalMean, level=0, period=(1992,2011), years_to_display=nothing)
default_options(::Val{:ECCO_map}) = (plot_type=:ECCO_map, statistic="mean", time=1)
default_options(::Val{:ECCO_DepthTime}) = (
    plot_type=:ECCO_DepthTime, period=(1992,2011), factor=1,
    level=1, klims=(1,50), years_to_display=nothing,
)
default_options(::Val{:ECCO_OHT1}) = (plot_type=:ECCO_OHT1, period=(1992,2011), years_to_display=nothing)
default_options(::Val{:ECCO_Overturn2}) = (plot_type=:ECCO_Overturn2, period=(1992,2011), years_to_display=nothing)
default_options(::Val{:ECCO_Overturn1}) = (plot_type=:ECCO_Overturn1, level=1, low1="auto", period=(1992,2011), years_to_display=nothing)
default_options(::Val{:ECCO_Transports}) = (plot_type=:ECCO_Transports, ncols=1, period=(1992,2011), years_to_display=nothing)

function finalize_options(o::NamedTuple)
    if haskey(o,:years_to_display) && isnothing(o.years_to_display) && haskey(o,:period)
        (y0,y1)=o.period
        o=merge(o,(years_to_display=(y0,y1+1),))
    end
    o
end

ECCOdiag(path::String,name::String,plot_type::Symbol; kwargs...) =
    ECCOdiag(path=path, name=name, options=setopt(plot_type;kwargs...))

"""
    setopt(plot_type::Symbol; kwargs...)

Construct a finalized options `NamedTuple` for `plot_type`, starting from
[`default_options`](@ref) and overriding/adding any fields passed as
keyword arguments.

```julia
setopt(plot_type::Symbol; kwargs...) =
    finalize_options(merge(default_options(plot_type), NamedTuple(kwargs)))
```

This is the preferred way to build an options object for an `ECCOdiag`/
`SSTdiag`, rather than writing `options=(...)` by hand — it guarantees
every plot-type-specific default is present, then applies
[`finalize_options`](@ref) to validate/normalize the result.

# Examples
```julia
o = setopt(:ECCO_OHT1; years_to_display=(1995,2005))
```

See also [`getopt`](@ref), [`default_options`](@ref), [`year_range`](@ref).
"""
setopt(plot_type::Symbol; kwargs...) =
    finalize_options(merge(default_options(plot_type), NamedTuple(kwargs)))

"""
    getopt(o::NamedTuple, k::Symbol, default)

Defensively read field `k` from the options `NamedTuple` `o`, returning
`default` instead of throwing a `KeyError` if `k` is absent.

```julia
getopt(o,k,default) = haskey(o,k) ? getproperty(o,k) : default
```

This is the preferred way to read optional fields from an `ECCOdiag`/`SSTdiag`
options object — in particular for fields (like `:years_to_display`) that may
not be present when the object was constructed by hand (e.g.
`options=(...)`) rather than via [`setopt`](@ref).

# Examples
```julia
getopt(o, :years_to_display, nothing)
```

See also [`setopt`](@ref), [`year_range`](@ref).
"""
getopt(o::NamedTuple,k::Symbol,default) = haskey(o,k) ? getproperty(o,k) : default

"""
    year_range(o::NamedTuple)

Return the `(Y0, Y1)` year range that should be used when averaging or
selecting a time sub-window for `o`, falling back to `o.period` when
`o.years_to_display` is absent or `nothing`.

```julia
(year0, year1) = o.period
yd = getopt(o, :years_to_display, nothing)
isnothing(yd) ? (year0, year1 + 1) : yd
```

This is the single defensive accessor for `years_to_display` across all
`ECCOdiag`/`SSTdiag` plot types, closing a `KeyError` risk for any options
object built by hand rather than via [`setopt`](@ref)/`finalize_options`.

Its meaning depends on the plot type's structure:
- For genuine time-series plot types (e.g. `ECCO_Overturn1`,
  `ECCO_Transports`, `ECCO_GlobalMean`, `ECCO_TimeLat`, `ECCO_DepthTime`),
  the returned `(Y0, Y1)` sets the plotted x-axis limits directly.
- For time-averaged profile/section plot types (`ECCO_OHT1`,
  `ECCO_Overturn2`), which plot an average against latitude/depth rather
  than time, `(Y0, Y1)` instead selects *which months get averaged* — the
  index offset is computed relative to `o.period`'s start year, not
  anchored at index 1:

```julia
(year0, year1) = o.period
(Y0, Y1) = year_range(o)
i0 = Int(round((Y0 - year0) * 12 + 1))
i1 = Int(round((Y1 - year0) * 12))
```

See also [`getopt`](@ref), [`setopt`](@ref), [`default_options`](@ref).
"""
function year_range(o::NamedTuple)
    (year0,year1)=o.period
    yd=getopt(o,:years_to_display,nothing)
    isnothing(yd) ? (year0,year1+1) : yd
end

import JLD2: load

load(x::ECCOdiag; file="",variable="single_stored_object") = begin
    if occursin("zonmean",x.name)
        fil=joinpath(x.path,x.name,"zonmean.jld2")
        fil=(ispath(fil) ? fil : joinpath(x.path,x.name,"zonmean2d.jld2"))
    elseif occursin("_glo2d",x.name)||occursin("_glo3d",x.name)
        fil=joinpath(x.path,x.name,"glo2d.jld2")
        fil=(ispath(fil) ? fil : joinpath(x.path,x.name,"glo3d.jld2"))
    elseif !isempty(file)
        fil=joinpath(x.path,x.name,file)
    else
        fil=joinpath(x.path,x.name,x.name*".jld2")
    end
    load(fil,variable)
end

export load

##

Base.@kwdef struct SSTdiag <: AbstractClimateDiagnostic
    path :: String = "unknown"
    name :: String = "unknown"
    options :: NamedTuple = NamedTuple()
    data :: AbstractArray = []
end

function SSTdiag(path::String, name::String, plot_type::Symbol; kwargs...)
    o = merge(default_options(Val(plot_type)), NamedTuple(kwargs))
    SSTdiag(path=path, name=name, options=o)
end

##

import DataFrames: DataFrame

Base.@kwdef struct SurfaceFluxDiag <: AbstractClimateDiagnostic
    options :: NamedTuple = NamedTuple()
    data :: Union{DataFrame,NamedTuple} = DataFrame()
end

##

Base.@kwdef struct SeaLevelAnomaly <: AbstractClimateDiagnostic
    path :: String = tempdir()
    name :: String = "unknown"
    options :: NamedTuple = NamedTuple()
    data :: AbstractArray = []
end

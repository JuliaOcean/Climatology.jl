
"""
    plot(x::SSTdiag)

Render the Makie figure appropriate for `x.options.plot_type`.

Dispatches on `string(x.options.plot_type)` to a method in `SST_plots`:

| `plot_type`         | renders via                        |
|:----------------------|:--------------------------------------|
| `"map_base"`           | `SST_plots.map_base()` (figure only)  |
| `"local_and_global"`   | `SST_plots.local_and_global(x)`       |
| `"by_year"`            | `SST_plots.by_year(x)`                |
| `"by_time"`            | `SST_plots.by_time(x)`                |
| `"TimeLat"`            | `SST_plots.TimeLat(x)`                |
| `"MHW"`                | `SST_plots.MHW(x)`                    |
| `"map"`                | `SST_plots.plot_sst_map(x)`           |
| `"ersst_map"`          | `SST_plots.ersst_map(x)`              |
| `"map_point"`          | `SST_plots.map_point(x)`              |

Throws an `ErrorException` if `x.options` is empty, or if `plot_type`
matches none of the cases above.
"""
	function plot(x::SSTdiag)
		isempty(x.options) && error("unknown options")
		o = x.options
		pt = string(o.plot_type)
		if pt=="map_base"
			fig,_,_ = SST_plots.map_base()
			fig
		elseif pt=="local_and_global"
			SST_plots.local_and_global(x)
		elseif pt=="by_year"
			SST_plots.by_year(x)
		elseif pt=="by_time"
			SST_plots.by_time(x)
		elseif pt=="TimeLat"
			SST_plots.TimeLat(x)
		elseif pt=="MHW"
			SST_plots.MHW(x)
		elseif pt=="map"
			SST_plots.plot_sst_map(x)
		elseif pt=="ersst_map"
			SST_plots.ersst_map(x)
		elseif pt=="map_point"
			SST_plots.map_point(x)
		else
			error("unknown plot_type")
		end
	end

module SST_plots

using Makie
import Climatology: load, Statistics, SSTdiag
import Climatology: MeshArrays, DataDeps
import Statistics: median

import Climatology: getopt
import Climatology.Dates: Date, year, dayofyear, daysinyear
import ClimatologyMakieExt: to_range!

# unified title resolution: explicit SSTdiag.name wins, else caller-supplied default
title_or(X::SSTdiag, default::AbstractString="") =
    (X.name in ("", "unknown") ? default : X.name)

#

"""
    SST_plots.decimal_years(ts, n)

Decimal year (e.g. `2024.5`) of each of the first `n` time steps of `ts`,
computed from `ts.year`, `ts.month` and `ts.day` (so it follows the actual
calendar, including leap years, however long the series is). Falls back to
a daily axis starting 1982-01-01 if `ts` has no date fields.
"""
function decimal_years(ts,n)
    if hasproperty(ts,:year) && hasproperty(ts,:month) && hasproperty(ts,:day)
        [y+(dayofyear(Date(y,m,d))-1)/daysinyear(y) for (y,m,d) in
            zip(ts.year[1:n],ts.month[1:n],ts.day[1:n])]
    else
        collect(1:n)/365.25 .+ 1982
    end
end

"""
    SST_plots.decimal_years(t::AbstractVector{<:Date})

Decimal year of each date in `t`: `year(t) + (dayofyear(t) - 1) / daysinyear(t)`,
so `2024-02-29` gives `2024 + 59/366`. Same values as the explorer's
`decimal_year`. The two-argument [`decimal_years`](@ref) is unchanged.
"""
decimal_years(t::AbstractVector{<:Date}) =
    [year(d) + (dayofyear(d)-1)/daysinyear(d) for d in t]

"""
    SST_plots.time_limits(o, tim)

x-axis limits `(year0, year1)` from option `period`; a `nothing` bound
(the default for `year1`) is taken from the data `tim` (decimal years), so
plots extend to the last available day.
"""
function time_limits(o,tim)
    y0,y1 = getopt(o,:period,(nothing,nothing))
    (isnothing(y0) ? floor(Int,tim[1]) : y0, isnothing(y1) ? tim[end] : y1)
end

"""
    SST_plots.year_ticks(year0, year1; step=5)

Tick positions at multiples of `step` years within `[year0, year1]`.
"""
year_ticks(year0,year1;step=5) = collect(ceil(Int,year0/step)*step : step : floor(Int,year1))

"""
    SST_plots.by_time(X::SSTdiag)

Plot the SST time series `X.options.timeseries`, optionally overlaid with
its seasonal climatology and/or anomaly.

Reads, via [`getopt`](@ref) (all optional, with the defaults shown):

- `X.options.timeseries` — a `NamedTuple` with fields `sst`, `clim`,
  `anom`, `title` (as produced by, e.g., `SST_timeseries.calc`)
- `show_anom::Bool = true` — overlay `timeseries.anom` in red
- `show_clim::Bool = true` — overlay `timeseries.clim` in orange
- `period::Tuple = (1982,nothing)` — x-axis limits; a `nothing` end
  (default) extends to the last available day

The time axis comes from `timeseries.year/month/day` (see
[`decimal_years`](@ref)); ticks are every 5 years.
"""
function by_time(X::SSTdiag)
    o = X.options
    ts = o.timeseries
    show_anom = getopt(o,:show_anom,true)
    show_clim = getopt(o,:show_clim,true)

    tim = decimal_years(ts,length(ts.sst))
    year0,year1 = time_limits(o,tim)
    f,a = lines(tim,ts.sst,label="SST",linewidth=4)
    show_clim ? lines!(a,tim,ts.clim,color=:orange,label="seasonal climatology",linewidth=1) : nothing
    show_anom ? lines!(a,tim,ts.anom,color=:red,label="SST - seasonal cycle") : nothing
    a.title = title_or(X, ts.title)
    a.xticks = year_ticks(year0,year1)
    xlims!(a,year0,year1)
    axislegend(a,position=:rb)
    f
end

"""
    SST_plots.by_year(X::SSTdiag)

Plot `X.options.timeseries.sst` as one line per calendar year against day
of year. Years are taken from `timeseries.year` (nothing is hardcoded), so the
highlighting follows the data: with `y1` the last year present,

- `y1` and `y1-1` (current and previous year): red, thicker
- `y1-4` to `y1-2` (the three years before that): blue
- all earlier years: gray
"""
function by_year(X::SSTdiag)
    ts = X.options.timeseries
    n = length(ts.sst)
    yr = ts.year[1:n]
    doy = [dayofyear(Date(y,m,d)) for (y,m,d) in zip(yr,ts.month[1:n],ts.day[1:n])]
    y1 = maximum(yr)

    f = Figure()
    a = Axis(f[1,1],xlabel="day of year")
    for (sel,col,lw) in ( (y->y<y1-4,:gray,1), (y->y1-4<=y<=y1-2,:blue,1), (y->y>=y1-1,:red,2) )
        for y in filter(sel,unique(yr))
            ii = findall(yr.==y)
            lines!(a,doy[ii],ts.sst[ii],color=col,linewidth=lw)
        end
    end
    a.title = title_or(X,"SST year by year (red=$(y1-1),$(y1); blue=$(y1-4):$(y1-2))")
    f
end

#

"""
    SST_plots.TimeLat(X::SSTdiag)

Plot a time-versus-latitude filled-contour anomaly diagram directly from
`X.options`.

Unlike the ECCO extension's [`TimeLat`](@ref) — which renders a
`NamedTuple` precomputed by a separate `ECCO_procs.TimeLat` step — this
method reads `X.options` directly, since SST currently has no
`_procs`-equivalent precompute stage. Fields read (via [`getopt`](@ref)
where optional):

- `X.options.timeseries` — used only for `.year`'s length, to build the
  time axis
- `X.options.zonal_mean` — the raw latitude × time array to contour
- `period::Tuple = (1982,nothing)` — x-axis limits (`nothing` end: last day of `zonal_mean`)
- `ylims::Tuple = (-90,90)` — y-axis (latitude) limits
- `clip_to_range::Bool = true` — clip the field to the fixed contour
  levels `(-2.0:0.25:2.0)/5` via `to_range!` before contouring

Latitude bins are inferred from `size(zonal_mean,1)`, assuming they evenly
tile `-90:90`. The time axis/field are subsampled by a factor of 7
(`x[1:7:end]`, `z[1:7:end,:]`) for plotting performance.
"""
function TimeLat(X::SSTdiag)
    o = X.options
    list = o.timeseries
    zm   = o.zonal_mean
    lat0,lat1   = getopt(o,:ylims,(-90,90))
    clip_to_range = getopt(o,:clip_to_range,true)

    x = decimal_years(list,size(zm,2))   # `list` may be longer than the zonal-mean table
    year0,year1 = time_limits(o,x)
    dy = Int(180/size(zm,1))
    y = collect(-90+dy/2:dy:90-dy/2)
    z = permutedims(zm)
    levs = (-2.0:0.25:2.0)/5.0

    clip_to_range ? to_range!(z,levs) : nothing
    fig1 = Figure(size = (900,400),markersize=0.1)
    ax1 = Axis(fig1[1,1], title=title_or(X,"OISST anomaly"),
        xticks=year_ticks(year0,year1),yticks=collect(-90.0:20.0:90.0),ylabel="latitude")
    hm1 = contourf!(ax1,x[1:7:end],y,z[1:7:end,:],levels=levs,colormap=:curl)
    Colorbar(fig1[1,2], hm1, height = Relative(0.65))
    xlims!(ax1,year0,year1)
    ylims!(ax1,lat0,lat1)
    fig1
end

#

function lowres_scatter(kdf,fig=[],ax=[]; input=[])
    (i,j) = ([x.i for x in kdf],[x.j for x in kdf])
    (ii,jj) = (10*i.-5,10*j.-95)
    if isa(fig,Array)
        f,a = scatter(ii,jj,color=input,markersize=10)
        c=(:blue,:red)
    else
        (f,a) = (ax,fig)
        c=(:skyblue,:pink)
    end
    text!(a,ii.+1,jj,text=string.(i),fontsize=11,color=c[1])
    text!(a,ii.+1,jj.-3,text=string.(j),fontsize=11,color=c[2])
    f
end

function local_and_global(X::SSTdiag)
    o = X.options
    ts        = o.timeseries
    ts_global = o.timeseries_global
    ylim0,ylim1 = getopt(o,:ylims,(-2.5,2.5))

    tim = decimal_years(ts,length(ts.anom))
    year0,year1 = time_limits(o,tim)
    fig,ax,li = lines(tim,ts.anom .-median(ts.anom),label="local")
    lines!(tim,ts_global.anom .-median(ts_global.anom),label="global")
    ax.title = title_or(X,"local and global SST anomalies")
    ax.xticks = year_ticks(year0,year1)
    xlims!(ax,year0,year1)
    ylims!(ax,ylim0,ylim1)
    axislegend(ax,position = :rb)
    fig
end

"""
    SST_plots.map_base()

Build a base `Figure`/`Axis` showing the Blue Marble Next Generation
basemap image, for use as a background under subsequent `heatmap!`/
`scatter!` calls (see [`plot_sst_map`](@ref)).

Returns `(fig, ax, im)`.

The image is fetched via `MeshArrays.mydatadep("basemap_jpg1")`, then
`reverse`d, `permutedims`d, and `circshift`ed by `(1800,0)` pixels to
align its native orientation and longitude origin with the `[0,360)`
convention used elsewhere in this package (displayed on `-0.05 .. 359.95`
in `x`, `-89.95 .. 89.95` in `y`). Axis decorations are hidden.
"""
function map_base()
    earth_jpg = joinpath(MeshArrays.mydatadep("basemap_jpg1"),
        "Blue_Marble_Next_Generation_+_topography_+_bathymetry.jpg")
    earth_img = load(earth_jpg)
    earth_img = reverse(permutedims(earth_img),dims=2)
    earth_img = circshift(earth_img,(1800,0))

    fig = with_theme(Figure,theme_light())
    ax = Axis(fig[1, 1])
    im = image!(ax, -0.05 .. 359.95, -89.95 .. 89.95, earth_img)
    hidedecorations!(ax)
    fig,ax,im
end

##

"""
    SST_plots.MHW(X::SSTdiag)

Plot the SST anomaly time series (`X.options.timeseries.sst .-
X.options.timeseries.clim`), highlighting marine heat wave (MHW) periods
in red.

A point is classified as an MHW day when the anomaly exceeds
`X.options.timeseries.high` — the day-of-year 90th percentile threshold
computed by `SST_timeseries.calc_quantile` (see that function for how
`.high` is derived); all other points are blue. `period`, read via
[`getopt`](@ref) with default `(1982,nothing)` (`nothing` end: last available
day), sets the x-axis limits.
"""
function MHW(X::SSTdiag)
    o = X.options
    ts = o.timeseries

    x = ts.sst-ts.clim
    y = fill(:blue,size(x))
    y[findall(x.>=ts.high)].=:red
    tim = decimal_years(ts,length(ts.sst))
    year0,year1 = time_limits(o,tim)

    fig,ax,li = lines(tim,x,color=y)
    ax.xticks = year_ticks(year0,year1)
    xlims!(ax,year0,year1)
    ax.title = title_or(X,"SST anomaly with extreme warm periods in red")
    fig
end

"""
    SST_plots.plot_sst_map(X::SSTdiag)

Plot an SST (or anomaly) map over the [`map_base`](@ref) basemap, with an
optional grid overlay and a highlighted point of interest.

Reads `X.options.map_data`, a `NamedTuple` expected to provide: `lon`,
`lat`, `field` (the gridded values to show via `heatmap!`), `colormap`,
`colorrange`, `showgrid::Bool` (whether to overlay [`lowres_scatter`](@ref)'s
coarse-grain grid-index labels), and `lon1`/`lat1` (coordinates of a
single point, marked with both a blue circle and a yellow X — e.g. to
indicate the location a companion time-series plot corresponds to).
"""
function plot_sst_map(X::SSTdiag)
    md = X.options.map_data
    fig,ax,_ = map_base()
    hm = heatmap!(ax,md.lon,md.lat,md.field,colormap=md.colormap,colorrange=md.colorrange)
    md.showgrid ? lowres_scatter(ax) : nothing
    scatter!(ax,md.lon1,md.lat1,marker=:circle,color=:blue,markersize=30)
    scatter!(ax,md.lon1,md.lat1,marker=:x,color=:yellow,markersize=15)
    Colorbar(fig[1, 2],hm)
    ax.title = title_or(X,"SST map")
    fig
end

"""
    SST_plots.ersst_map(X::SSTdiag)

Plot an ERSST monthly anomaly map over the [`map_base`](@ref) basemap, from
`X.options.map_data` as returned by [`SST_FILES.ersst_map_data`](@ref). Unlike
[`plot_sst_map`](@ref), it draws no point marker and no grid overlay, so it
needs only `field`, `lon`, `lat`, `colorrange`, `colormap` and `title`.
The title is `X.name` when set, else `map_data.title`.
"""
function ersst_map(X::SSTdiag)
    md = X.options.map_data
    fig,ax,_ = map_base()
    hm = heatmap!(ax,md.lon,md.lat,md.field,colormap=md.colormap,colorrange=md.colorrange)
    ax.title = title_or(X,md.title)
    Colorbar(fig[1, 2],hm)
    fig
end

"""
    SST_plots.map_point(X::SSTdiag)

Map from `X.options.map_data` (as in [`plot_sst_map`](@ref), but with its own
colour settings and no grid overlay) next to the 10° box series for a point.

Options: `point` (from [`SST_coarse_grain.point_series`](@ref), or `nothing` over
land), `timeseries` (from `SST_timeseries.calc`, or `nothing`), `lon`, `lat`
(requested point, shown as a black cross), `show_anom` (plot `sst - clim` instead
of `sst`). Over land the right panel says so instead of drawing a series.
The left title is `X.name` when set, else `map_data.title`.
"""
function map_point(X::SSTdiag)
    o = X.options
    m = o.map_data
    ps = getopt(o,:point,nothing)
    ts = getopt(o,:timeseries,nothing)
    lon_req = getopt(o,:lon,205.0)
    lat_req = getopt(o,:lat,25.0)
    show_anom = getopt(o,:show_anom,false)

    fig = Figure(size = (1300, 480))
    ax1 = Axis(fig[1, 1]; title = title_or(X,m.title), xlabel = "longitude (°E)",
               ylabel = "latitude (°N)", aspect = DataAspect())
    z = [ismissing(v) ? NaN32 : Float32(v) for v in m.field]
    hm = heatmap!(ax1, m.lon, m.lat, z; colorrange = m.colorrange, colormap = m.colormap,
                  nan_color = :gray80)
    Colorbar(fig[1, 2], hm; label = m.metadata.units)
    limits!(ax1, 0, 360, -90, 90)
    scatter!(ax1, [Float64(lon_req)], [Float64(lat_req)]; marker = :xcross, color = :black, markersize = 14)

    ax2 = Axis(fig[1, 3]; xlabel = "year")
    if ps === nothing
        ax2.title = "No ocean cell at ($(lon_req)°E, $(lat_req)°N)"
        hidedecorations!(ax2)
    else
        h = ps.metadata.dlon / 2
        lines!(ax1, ps.lon .+ h .* [-1, 1, 1, -1, -1], ps.lat .+ h .* [-1, -1, 1, 1, -1];
               color = :black, linewidth = 2)
        x = decimal_years(ps.time)
        if show_anom && ts !== nothing && length(ts.sst) == length(x)
            lines!(ax2, x, ts.sst .- ts.clim; linewidth = 0.5)
            hlines!(ax2, 0; color = :gray50)
            ax2.ylabel = "SST anomaly (degC, vs 1992–2011)"
        else
            lines!(ax2, x, ps.sst; linewidth = 0.5)
            ax2.ylabel = "SST (degC)"
        end
        ax2.title = "$(ps.metadata.dlon)°×$(ps.metadata.dlon)° box mean, centre ($(ps.lon)°E, $(ps.lat)°N)"
    end
    colsize!(fig.layout, 3, Relative(0.4))
    fig
end

end

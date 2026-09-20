
import Climatology: year_range

"""
    plot(x::ECCOdiag)

Render the Makie figure appropriate for `x.options.plot_type`.

Dispatches on `x.options.plot_type` (a `Symbol`, see [`default_options`](@ref)),
combined for a few plot types with `x.name`:

| `plot_type`          | requires `x.name ==` | renders via                                              |
|:----------------------|:-----------------------|:-------------------------------------------------------------|
| `:ECCO_map`            |                         | [`ECCO_map`](@ref) applied to `ECCO_procs.ECCO_map(x)`        |
| `:ECCO_TimeLat`        |                         | [`TimeLat`](@ref) applied to `ECCO_procs.TimeLat(x)`          |
| `:ECCO_TimeLatAnom`    |                         | [`TimeLat`](@ref) applied to `ECCO_procs.TimeLat(x)`          |
| `:ECCO_DepthTime`      |                         | [`DepthTime`](@ref) applied to `ECCO_procs.DepthTime(x)`      |
| `:ECCO_GlobalMean`     |                         | [`glo`](@ref) applied to `ECCO_procs.glo(x)`                  |
| `:ECCO_OHT1`           | `"OHT"`                 | [`OHT`](@ref)                                                 |
| `:ECCO_Overturn2`      | `"overturn"`            | [`OverturnStreamfunction`](@ref) (formerly `figov2`)          |
| `:ECCO_Overturn1`      | `"overturn"`            | [`OverturnTimeseries`](@ref) (formerly `figov1`)              |
| `:ECCO_Transports`     | `"trsp"`                | [`transport`](@ref)                                           |

For the first five rows, `x` (an `ECCOdiag`) is first passed through the
same-named function in `ECCO_procs` to compute a plain `NamedTuple` of plot
data, which is then rendered by the like-named method in this extension —
i.e. `ECCO_map`/`TimeLat`/`DepthTime`/`glo` each name *two* functions
(one data-producing, in `ECCO_procs`; one figure-producing, here). The
remaining plot types (`OHT`, `OverturnTimeseries`, `OverturnStreamfunction`,
`transport`) instead read and process data directly from `x.path`/`x.name`,
without a separate `ECCO_procs` precompute step.

Prints `"unknown option (a)"` if `x.options` is empty, or
`"unknown option (b)"` if `plot_type`/`name` match none of the cases
above, rather than throwing an error.
"""
function plot(x::ECCOdiag)
	if !isempty(x.options)
		o=x.options
		pt=string(o.plot_type)
		if pt=="ECCO_map"
			ECCO_map(ECCO_procs.ECCO_map(x))			
		elseif pt in ("ECCO_TimeLat","ECCO_TimeLatAnom")
    		TimeLat(ECCO_procs.TimeLat(x))
		elseif pt=="ECCO_DepthTime"
			DepthTime(ECCO_procs.DepthTime(x))
		elseif pt=="ECCO_GlobalMean"
			glo(ECCO_procs.glo(x))
		elseif x.name=="OHT" && pt=="ECCO_OHT1"
			OHT(x)
		elseif x.name=="overturn" && pt=="ECCO_Overturn2"
			OverturnStreamfunction(x)
		elseif x.name=="overturn" && pt=="ECCO_Overturn1"
			OverturnTimeseries(x)
		elseif x.name=="trsp" && pt=="ECCO_Transports"
			transport(x)
		else
			println("unknown option (b)")	
		end
	else
		println("unknown option (a)")
	end
end

##

#	years_to_display=(1960,2023)
years_to_display=(1980,2024)

"""
    time_average_indices(o::NamedTuple, nt::Int)

Compute the clamped `(i0,i1)` month-index range to average over, given
options `o` (using `o.period` and [`year_range`](@ref)) and the actual
number of time records `nt` available in the loaded data.

```julia
(year0,year1) = o.period
(Y0,Y1) = year_range(o)
i0 = Int(round((Y0-year0)*12+1))
i1 = Int(round((Y1-year0)*12))
```

`i0` and `i1` are clamped to `1:nt` (rather than left to error via
`BoundsError`) so that a `years_to_display`/`period` window extending
past the end (or before the start) of the loaded record silently
truncates to the available data instead of failing. A warning is emitted
via `@warn` when clamping actually occurs, to surface a likely
mismatched-window mistake without hard-failing the plot.
"""
function time_average_indices(o::NamedTuple, nt::Int)
    (year0,year1) = o.period
    (Y0,Y1) = year_range(o)
    i0 = Int(round((Y0-year0)*12+1))
    i1 = Int(round((Y1-year0)*12))
    if i0 < 1 || i1 > nt
        @warn "requested averaging window ($Y0,$Y1) exceeds available data range; clamping"
    end
    (max(i0,1), min(i1,nt))
end

##

function axtr1(ax,namtr,pth_out,list_trsp,year0,year1;years_to_display=years_to_display)
	itr=findall(list_trsp.==namtr)[1]
	tmp=vec(load(ECCOdiag(path=pth_out,name="trsp")))[itr]
	
	nt=size(tmp.val,2)
	x=vec(0.5:nt)

	txt=tmp.nam[1:end-5]
	val=1e-6*vec(sum(tmp.val,dims=1)[:])
	valsmo = runmean(val, 12)

	x=vec(0.5:nt)
	x=year0 .+ x./12.0

	hm1=lines!(ax,x,val,label="ECCO estimate")
	valsmo[1:5].=NaN
	valsmo[end-4:end].=NaN
	lines!(ax,x,valsmo,linewidth=4.0,color=:red)
	xlims!(ax,years_to_display)
end

"""
    transport(X::ECCOdiag)

Plot 12-month running-mean volume transport (in Sv) time series, one
panel per named section in `X.options.namtrs`.

Reads from `X.options`: `namtrs` (section names to plot, matched against
`list_trsp`), `ncols` (panels per row), `list_trsp` (full list of section
names as stored in the `"trsp"` data file), and `period` (a `(year0,year1)`
tuple used for tick spacing). The plotted x-axis range is
`years_to_display = year_range(X.options)`.

Each panel is rendered by the internal helper `axtr1`, which loads and
converts the corresponding section's transport from `X.path`.
"""
function transport(X::ECCOdiag)
    o=X.options
    namtrs=o.namtrs
    ncols=o.ncols
    list_trsp=o.list_trsp
    (year0,year1)=o.period
    years_to_display=year_range(o)
    pth_out=X.path

    fig1 = ncols>1 ? Figure(size=(2000,1000),markersize=0.1) : Figure(size=(900,400),markersize=0.1)
    for na in 1:length(namtrs)
        txt=namtrs[na]
        jj=div.(na,ncols,RoundUp)
        kk=na-(jj.-1)*ncols
        ax1 = Axis(fig1[jj,kk], title=" $txt (in Sv)",
            xticks=(year0:4:year1),ylabel="transport, in Sv")
        axtr1(ax1,namtrs[na],pth_out,list_trsp,year0,year1,years_to_display=years_to_display)
    end
    fig1
end

"""
    OverturnTimeseries(X::ECCOdiag)

Plot 12-month running-mean overturning transport (in Sv) time series at a
fixed depth level, for a fixed set of latitudes.

`X.options.plot_type == :ECCO_Overturn1`, `X.name == "overturn"`. Loads
`X.name` data from `X.path`, then plots one line per latitude index in
`115:10:145` against `lats = -89.0:89.0` — i.e. **25°N, 35°N, 45°N,
55°N** (this index range must be kept in sync with `lats` if either is
changed). `X.options.level` selects the depth index (shown as `kk` in the
title); `X.options.low1` overrides the y-axis lower bound (`"auto"` uses
the fixed `(5,20)` Sv range).

The x-axis is set to `years_to_display = year_range(X.options)`, which
for this plot type is a genuine axis range (contrast with [`OHT`](@ref)
and [`OverturnStreamfunction`](@ref), where the same accessor instead selects a
time-averaging window).
"""
function OverturnTimeseries(X::ECCOdiag)
    o=X.options
    level=o.level
    low1=o.low1
    (year0,year1)=o.period
    years_to_display=year_range(o)

    tmp=-1e-6*load(ECCOdiag(path=X.path,name=X.name))
    nt=size(tmp,3)
    x=vec(0.5:nt)
    x=year0 .+ x./12.0
    lats=vec(-89.0:89.0)

    fig1 = Figure(size = (900,400),markersize=0.1)
    ax1 = Axis(fig1[1,1],ylabel="Sv",
        title="Global Overturning, in Sv, at kk=$(level)",
        xticks=(year0:4:year1))
    for ll in 115:10:145
        ov=tmp[ll,level,:]
        ov=runmean(ov, 12)
        ov[1:5].=NaN
        ov[end-4:end].=NaN
        lines!(x,ov,label="$(lats[ll])N")
    end
    xlims!(ax1,years_to_display)
    ylims!(ax1,(5,20))
    low1!="auto" ? ylims!(ax1,(low1,20.0)) : nothing
    fig1[1, 2] = Legend(fig1, ax1, "estimate", framevisible = false)
    fig1
end

"""
    OverturnStreamfunction(X::ECCOdiag; ClipToRange=true)

Plot the time-averaged meridional overturning streamfunction (in Sv) as a
filled contour over latitude and depth.

Like [`OHT`](@ref), and unlike the `ECCO_procs`-mediated plot types,
`OverturnStreamfunction` takes `X` directly (`X.options.plot_type == :ECCO_Overturn2`,
`X.name == "overturn"`), loading `X.name` data from `X.path` and
computing the time average inline over month indices `i0:i1` derived from
`year_range(X.options)` relative to `X.options.period` (see [`OHT`](@ref)
for the exact index formula, and [`year_range`](@ref) for why this plot
type's `years_to_display` selects an averaging window rather than an
axis range).

`X.options.grid` supplies the vertical grid (`Γ.RF`) used for the depth
axis.

When `ClipToRange` is `true` (default), the time-averaged field is
clipped to the fixed contour levels `-40:5:40` Sv via `to_range!` before
contouring.
"""
function OverturnStreamfunction(X::ECCOdiag; ClipToRange=true)
    o=X.options
    Γ=o.grid
    (year0,year1)=o.period
    (Y0,Y1)=year_range(o)

    tmp=-1e-6*load(ECCOdiag(path=X.path,name=X.name))
    i0,i1=time_average_indices(o,size(tmp,3))
    ovmean=dropdims(mean(tmp[:,:,i0:i1],dims=3),dims=3)
    x=vec(-89.0:89.0); y=reverse(vec(Γ.RF[1:end-1]))
    z=reverse(ovmean,dims=2); z[z.==0.0].=NaN
    levs=(-40.0:5.0:40.0)
    ClipToRange ? to_range!(z,levs) : nothing
    fig1 = Figure(size = (900,400),markersize=0.1)
    ax1 = Axis(fig1[1,1], title="Meridional Overturning Streamfunction (in Sv, $(Y0)-$(Y1-1) mean)",
            xlabel="latitude",ylabel="depth (in m)")
    hm1=contourf!(ax1,x,y,z,levels=levs)
    Colorbar(fig1[1,2], hm1, height = Relative(0.65))
    fig1
end

"""
    OHT(X::ECCOdiag)

Plot the time-averaged northward ocean heat transport (in PW) versus
latitude.

Unlike [`ECCO_map`](@ref)/[`TimeLat`](@ref)/[`DepthTime`](@ref)/[`glo`](@ref),
`OHT` takes `X` (`X.options.plot_type == :ECCO_OHT1`, `X.name == "OHT"`)
directly rather than a precomputed `NamedTuple`: it loads `"MHT"` data
from `X.path` and computes the plotted time average inline.

The averaging window is `(Y0,Y1) = year_range(X.options)`, converted to
month indices relative to `X.options.period`'s start year:

```julia
(year0,year1) = X.options.period
(Y0,Y1) = year_range(X.options)
i0 = Int(round((Y0-year0)*12+1))
i1 = Int(round((Y1-year0)*12))
```

For this plot type, `years_to_display`/[`year_range`](@ref) selects
*which months get averaged*, not a plotted x-axis range — the x-axis here
is latitude, not time.
"""
function OHT(X::ECCOdiag)
    o=X.options
    (year0,year1)=o.period
    (Y0,Y1)=year_range(o)
    pth_out=X.path

    tmp=load(ECCOdiag(path=pth_out,name="MHT"))
    i0,i1=time_average_indices(o,size(tmp,2))
    MT=vec(mean(tmp[:,i0:i1],dims=2))

    x=vec(-89.0:89.0)
    fig1 = Figure(size = (900,400),markersize=0.1)
    ax1 = Axis(fig1[1,1], title="Northward Heat Transport (in PW, $(Y0)-$(Y1-1) mean)",
        xticks=(-90.0:10.0:90.0),yticks=(-2.0:0.25:2.0),
        xlabel="latitude",ylabel="Transport (in PW)")
    lines!(x,MT)
    ylims!(ax1,(-2.0,2.0))
    fig1
end

"""
code snippet for conversion to ZJoule from older version of `glo(gl1)`

```
	if false
		fac=4e6*1.335*10^9*10^9/1e21
		ttl="Ocean Heat Uptake (Zetta-Joules)"
		zlb="Zetta-Joules"
		rng=(-100.0,300.0)
		y=fac*(gl1.y.-gl1.y[1])
	else
		y=gl1.y
	end
```
"""

"""
    glo(gl1)

Render a global-mean time series line plot from `gl1`.

`gl1` is the `NamedTuple` returned by `ECCO_procs.glo` — not an
`ECCOdiag` directly (see [`plot`](@ref)) — with fields `x` (time), `y`
(global-mean value), `txt` (used as both title and y-axis label), `rng`
(y-axis limits), `year0`/`year1` (x-tick spacing), and
`years_to_display` (x-axis limits).
"""
function glo(gl1)
    fig1 = Figure(size = (900,400),markersize=0.1)
    ax1 = Axis(fig1[1,1], title="Global Mean $(gl1.txt)",
        xticks=collect(gl1.year0:4:gl1.year1),ylabel=gl1.txt)
    lines!(ax1,gl1.x,gl1.y)
    xlims!(ax1,gl1.years_to_display)
    ylims!(ax1,gl1.rng)
    fig1
end

"""
    DepthTime(XYZ; ClipToRange=true)

Render a time-versus-depth filled-contour diagram from `XYZ`.

`XYZ` is the `NamedTuple` returned by `ECCO_procs.DepthTime` — not an
`ECCOdiag` directly (see [`plot`](@ref)) — with fields `x` (time), `y`
(depth), `z` (depth × time field), `levels`, `title`, `ylims` (a
`(depth0,depth1)` tuple, typically given reversed so depth increases
downward on screen), and `years_to_display` (x-axis limits — see
[`year_range`](@ref)).

When `ClipToRange` is `true` (default), `XYZ.z` is clipped in place to
`XYZ.levels`'s range via `to_range!` before contouring.
"""
function DepthTime(XYZ; ClipToRange=true)
    ClipToRange ? to_range!(XYZ.z,XYZ.levels) : nothing
    fig1 = Figure(size=(900,400),markersize=0.1)
    ax1 = Axis(fig1[1,1], title=XYZ.title,
        xticks=collect(XYZ.year0:4:XYZ.year1))
    hm1=contourf!(ax1,XYZ.x,XYZ.y,XYZ.z,levels=XYZ.levels,colormap=:turbo)
    Colorbar(fig1[1,2], hm1, height=Relative(0.65))
    xlims!(ax1,XYZ.years_to_display)
    ylims!(ax1,XYZ.ylims)
    fig1
end

"""
    TimeLat(XYZ; ClipToRange=true)

Render a time-versus-latitude filled-contour (Hovmöller) diagram from
`XYZ`.

`XYZ` is the `NamedTuple` returned by `ECCO_procs.TimeLat` — not an
`ECCOdiag` directly (see [`plot`](@ref)) — with fields `x` (time), `y`
(latitude), `z` (time × latitude field), `levels`, `title`, `ylims` (a
`(lat0,lat1)` tuple), and `years_to_display` (an `(x0,x1)` tuple used for
the x-axis limits — see [`year_range`](@ref)).

When `ClipToRange` is `true` (default), `XYZ.z` is clipped in place to
`XYZ.levels`'s range via `to_range!` before contouring.
"""
function TimeLat(XYZ; ClipToRange=true)
    ClipToRange ? to_range!(XYZ.z,XYZ.levels) : nothing
    fig1 = Figure(size=(900,400),markersize=0.1)
    ax1 = Axis(fig1[1,1], title=XYZ.title,
        xticks=collect(XYZ.year0:4:XYZ.year1),yticks=collect(-90.0:20.0:90.0),ylabel="latitude")
    hm1=contourf!(ax1,XYZ.x,XYZ.y,XYZ.z,levels=XYZ.levels,colormap=:turbo)
    Colorbar(fig1[1,2], hm1, height=Relative(0.65))
    xlims!(ax1,XYZ.years_to_display)
    ylims!(ax1,XYZ.ylims...)
    fig1
end

"""
    ECCO_map(X; ClipToRange=true)

Render a filled-contour map of `X.field` over the longitude/latitude grid
`X.λ`.

`X` is the `NamedTuple` returned by `ECCO_procs.ECCO_map` — not an
`ECCOdiag` directly (see [`plot`](@ref) for how the two connect) — with
fields `λ` (interpolation target, providing `lon`/`lat`), `field`,
`levels`, and `title`.

When `ClipToRange` is `true` (default), values in `X.field` outside the
range of `X.levels` are clipped in place via `to_range!` before
contouring, avoiding blank contour gaps for out-of-range data.
"""
function ECCO_map(X; ClipToRange=true)
	ClipToRange ? to_range!(X.field,X.levels) : nothing
	fig = Figure(size = (900,600), backgroundcolor = :grey95)
	ax = Axis(fig[1,1], title=X.title,xlabel="longitude",ylabel="latitude")
	hm1=contourf!(ax,X.λ.lon[:,1],X.λ.lat[1,:],X.field,levels=X.levels,colormap=:turbo)
	Colorbar(fig[1,2], hm1, height = Relative(0.65))
	fig	
end

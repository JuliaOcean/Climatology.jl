module ClimatologyMakieExt
	using Makie, Climatology
	import Climatology: Statistics, RollingFunctions, plot_examples, load
	import Climatology: ECCOdiag, SSTdiag, SeaLevelAnomaly, SurfaceFluxDiag

	import Statistics: mean
	import Makie: plot
	import RollingFunctions: runmean

"""
    to_range!(DD, levs)
    to_range!(DD, levs::Tuple)

Clip array `DD` in place so that all values fall strictly within the
range spanned by contour levels `levs`, avoiding blank/uncontoured
regions where `DD` lies exactly at or beyond `levs`'s extremes.

Values `<= levs[1]` are set to `levs[1] + (levs[2]-levs[1])/100`; values
`>= levs[end]` are set to `levs[end] - (levs[end]-levs[end-1])/100` — i.e.
nudged just inside the boundary rather than removed or `NaN`-filled.

The `Tuple` method accepts a plain `(lo,hi)` range (as used for, e.g.,
`ECCO_map`'s `figov2` fixed `-40:5:40` levels) and expands it to 10
evenly-spaced levels via `range(levs[1],levs[2],length=10)` before
delegating to the vector method.

Used by [`ECCO_map`](@ref), [`TimeLat`](@ref), [`DepthTime`](@ref),
[`figov2`](@ref) (ECCO), and `SST_plots.TimeLat` (SST) prior to
`contourf!`.
"""
	function to_range!(DD,levs)
		DD[findall(DD.<=levs[1])].=levs[1]+(levs[2]-levs[1])/100
		DD[findall(DD.>=levs[end])].=levs[end]-(levs[end]-levs[end-1])/100
	end

	to_range!(DD,levs::Tuple) = to_range!(DD,range(levs[1],levs[2],length=10))

	include("Makie/ECCO.jl")
	include("Makie/OISST.jl")
	include("Makie/SLA.jl")
	include("Makie/SurfaceFluxes.jl")
end
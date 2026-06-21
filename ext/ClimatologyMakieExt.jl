module ClimatologyMakieExt
	using Makie, Climatology
	import Climatology: Statistics, RollingFunctions, plot_examples, load
	import Climatology: ECCOdiag, SSTdiag, SeaLevelAnomaly, SurfaceFluxDiag

	import Statistics: mean
	import Makie: plot
	import RollingFunctions: runmean

	include("Makie/ECCO.jl")
	include("Makie/OISST.jl")
	include("Makie/SLA.jl")
	include("Makie/SurfaceFluxes.jl")
end

##

abstract type AbstractClimateDiagnostic <: Any end

default_options(plot_type::Symbol) = default_options(Val(plot_type))
default_options(::Val{T}) where T = (plot_type=T,)

##

Base.@kwdef struct ECCOdiag <: AbstractClimateDiagnostic
    path :: String = tempdir()
    name :: String = "unknown"
    options :: NamedTuple = NamedTuple()
    data :: AbstractArray = []
end

###

default_options(::Val{T}) where T = (plot_type=T, period=(1982,2024))

default_options(::Val{:by_year})          = (plot_type=:by_year,)
default_options(::Val{:by_time})          = (plot_type=:by_time, period=(1982,2024), show_anom=true, show_clim=true)
default_options(::Val{:MHW})              = (plot_type=:MHW, period=(1982,2024))
default_options(::Val{:local_and_global}) = (plot_type=:local_and_global, period=(1982,2024), ylims=(-2.5,2.5))
default_options(::Val{:map_base})         = (plot_type=:map_base,)
default_options(::Val{:map})              = (plot_type=:map,)
default_options(::Val{:TimeLat})          = (plot_type=:TimeLat, period=(1982,2024), ylims=(-90,90), clip_to_range=true)

default_options(::Val{:ECCO_TimeLat}) = (
    plot_type=:ECCO_TimeLat, select_method=0, period=(1992,2011),
    level=1, ylims=(-90,90), colormap_factor=1, years_to_display=nothing,
)
default_options(::Val{:ECCO_TimeLatAnom}) =
    merge(default_options(Val(:ECCO_TimeLat)), (plot_type=:ECCO_TimeLatAnom, select_method=1))
default_options(::Val{:ECCO_Overturn2}) = (plot_type=:ECCO_Overturn2,)
default_options(::Val{:ECCO_GlobalMean}) = (plot_type=:ECCO_GlobalMean, level=0, period=(1992,2011), years_to_display=nothing)
default_options(::Val{:ECCO_map}) = (plot_type=:ECCO_map, statistic="mean", time=1)
default_options(::Val{:ECCO_DepthTime}) = (
    plot_type=:ECCO_DepthTime, period=(1992,2011), factor=1,
    level=1, klims=(1,50), years_to_display=nothing,
)
default_options(::Val{:ECCO_OHT1}) = (plot_type=:ECCO_OHT1,)
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

setopt(plot_type::Symbol; kwargs...) =
    finalize_options(merge(default_options(plot_type), NamedTuple(kwargs)))

getopt(o::NamedTuple,k::Symbol,default) = haskey(o,k) ? getproperty(o,k) : default

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


##

abstract type AbstractClimateDiagnostic <: Any end

##

Base.@kwdef struct ECCOdiag <: AbstractClimateDiagnostic
    path :: String = tempdir()
    name :: String = "unknown"
    options :: NamedTuple = NamedTuple()
    data :: AbstractArray = []
end

ECCOdiag(path::String,name::String,plot_type::Symbol; kwargs...) =
    ECCOdiag(path=path, name=name, options=setopt(plot_type;kwargs...))

setopt(plot_type::Symbol; kwargs...) =
    ECCO_procs.finalize_options(merge(ECCO_procs.default_options(plot_type), NamedTuple(kwargs)))

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

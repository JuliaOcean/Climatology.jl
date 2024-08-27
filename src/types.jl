
##

abstract type AbstractClimateDiagnostic <: Any end

struct ECCOdiag <: AbstractClimateDiagnostic
    path::String
    name::String
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



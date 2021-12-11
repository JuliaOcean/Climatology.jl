
using MeshArrays, MITgcmTools, NCTiles
using JLD2, UUIDs, Unitful, Printf, TOML
using Distributed, SharedArrays

include("ECCO_helper_functions.jl")

γ,Γ=GridLoad_Main()
nr=length(Γ.DRF)
LC=LatitudeCircles(-89.0:89.0,Γ)
nl=length(LC)

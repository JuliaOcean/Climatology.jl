
using OceanStateEstimation, Printf, TOML, MeshArrays

γ,Γ,LC=ECCO_helper_functions.GridLoad_Main()
nr=length(Γ.DRF)
nl=length(LC)

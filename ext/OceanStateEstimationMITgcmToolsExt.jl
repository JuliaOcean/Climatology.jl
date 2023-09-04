module OceanStateEstimationMITgcmToolsExt

    import OceanStateEstimation: read_nctiles_alias, untargz_alias 
    
    import MITgcmTools: read_nctiles, untargz

    read_nctiles_alias(args...;kwargs...)=read_nctiles(args...;kwargs...)

    untargz_alias(fil)=untargz(fil)
    
end



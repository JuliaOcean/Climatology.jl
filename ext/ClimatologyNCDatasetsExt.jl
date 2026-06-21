module ClimatologyNCDatasetsExt

    import Climatology: ECCO, load, read_Dataset, ECCOdiags_to_nc
    import Climatology: write_SST_climatology, SST_demo_path, to_monthly_file
    import Climatology: write_SLA_PODAAC, write_SLA_CMEMS
    import MeshArrays, Printf 
    import MeshArrays: GridSpec, Tiles, GridLoadVar    
    import NCDatasets: Dataset, defDim, defVar

    read_Dataset(args...;kwargs...)=Dataset.(args...;kwargs...)

    include("NCDatasets/ECCO.jl")
    include("NCDatasets/OISST.jl")
    include("NCDatasets/OHC.jl")
    include("NCDatasets/SLA.jl")

end

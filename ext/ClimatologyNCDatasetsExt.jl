module ClimatologyNCDatasetsExt

    import Climatology: read_Dataset 
    
    import NCDatasets: Dataset

    read_Dataset(args...;kwargs...)=Dataset.(args...;kwargs...)
    
end



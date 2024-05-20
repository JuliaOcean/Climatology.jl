module OceanStateEstimation

pkg_pth=dirname(pathof(OceanStateEstimation))

#read_Dataset : Placeholder to allow NCDatasets extension, which is activated by `using NCDatasets`.
#read_nctiles_alias : Placeholder to allow MITgcmTools extension, which is activated by `using MITgcmTools`.

function read_Dataset end
function read_nctiles_alias end
function read_mdsio_alias end

include("downloads.jl")
include("ECCO.jl")

import OceanStateEstimation.downloads: get_ecco_files, get_ecco_variable_if_needed, get_ecco_velocity_if_needed
import OceanStateEstimation.downloads: get_occa_variable_if_needed, get_occa_velocity_if_needed
import OceanStateEstimation.downloads: ECCOdiags_add, CBIOMESclim_download, MITPROFclim_download
import DataDeps; import DataDeps: @datadep_str

using Glob

"""
   examples()

List of examples provided in OceanStateEstimation.jl (full paths)
"""
function examples()
    nb=joinpath(abspath("/"),split(pathof(OceanStateEstimation),"/")[2:end-2]...,"examples")
#    ex=glob("*/*.jl",nb)
    ex_known=("CBIOMES_climatology_plot.jl","ECCO_standard_plots.jl",
    "HadIOD_viz.jl","NSLCT_notebook.jl","OptimalTransport_demo.jl")
    ex=[glob("*/"*e,nb)[1] for e in ex_known]
end

export @datadep_str, ECCOdiags_add
export get_ecco_variable_if_needed, get_ecco_velocity_if_needed
export get_occa_variable_if_needed, get_occa_velocity_if_needed

export ECCO, ECCO_helpers, ECCO_io, ECCO_diagnostics
export ScratchSpaces

__init__() = begin
    ScratchSpaces.__init__scratch()
    downloads.__init__standard_diags()
end

end # module
